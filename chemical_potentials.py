import sys
from itertools import combinations
from typing import IO, Tuple
from dataclasses import dataclass
from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import ArrayLike, NDArray
from numpy.lib.stride_tricks import sliding_window_view


@dataclass
class Filter(ABC):

    @abstractmethod
    def employ(self, data: ArrayLike) -> Tuple[ArrayLike, float]:
        pass
        
        
@dataclass
class IdentityFilter(Filter):

    def employ(self, data: ArrayLike) -> Tuple[ArrayLike, float]:
        return data, 1.0
        
        
@dataclass
class OneDimChunkFilter(Filter):

    axis: int
    low: float
    high: float

    def employ(self, data: ArrayLike) -> ArrayLike:
    
        mask = np.logical_and(self.low <= data[:, -3 + self.axis], data[:, -3 + self.axis] <= self.high)
        new_data = data[mask]
        
        return new_data, new_data.shape[0] / data.shape[0]


def get_chemical_potentials(data: ArrayLike, filter_: callable = IdentityFilter()) -> NDArray:

    data, p = filter_.employ(data)

    num_types = data.shape[1] - 4
    
    init_types = data[:, 0].astype(int)
    occupying_energies = data[:, 1:num_types + 1]
    occupying_energies *= p
    
    types = np.arange(num_types, dtype=int)
    type_pairs = list(combinations(types, r=2))
    num_pairs = len(type_pairs)
    
    coefficient_matrix = np.zeros((num_pairs + 1, num_types))
    b = np.zeros(num_pairs + 1)
    
    for index, (t1, t2) in enumerate(type_pairs):
        coefficient_matrix[index, t1] = 1.0
        coefficient_matrix[index, t2] = -1.0
        b[index] = np.mean(occupying_energies[:, t1] - occupying_energies[:, t2])
        
    coefficient_matrix[num_pairs, :] = np.mean(init_types - 1 == types[:, None], axis=1)
    
    reference_occupying_energies = occupying_energies[np.arange(occupying_energies.shape[0]), init_types - 1]
    b[num_pairs] = np.mean(reference_occupying_energies) / occupying_energies.shape[0]
    
    return np.linalg.lstsq(coefficient_matrix, b, rcond=None)[0]


def main():

    axes = {
        'x': 0,
        'y': 1,
        'z': 2
    }

    data = np.loadtxt(sys.argv[1])
    data = data[data[:, -1].argsort()]
    z_coords = data[:, -3 + axes[sys.argv[2]]]
    num_atoms_per_chunk = int(sys.argv[3])
    windows = sliding_window_view(z_coords, window_shape=num_atoms_per_chunk)
    window_centers = np.zeros(windows.shape[0])
    chemical_potentials = np.zeros((windows.shape[0], data.shape[1] - 4))
    for i, w in enumerate(windows):
        window_centers[i] = np.mean(w)
        filter_ = OneDimChunkFilter(axis=2, low=np.min(w), high=np.max(w))
        chemical_potentials[i, :] = get_chemical_potentials(data, filter_=filter_)
        
    header='center'
    for i in np.arange(data.shape[1] - 4, dtype=int):
        header += f' mu_{i + 1:.0f}'
        
    np.savetxt(
        sys.argv[4],
        np.vstack((window_centers, *chemical_potentials.T)).T,
        header=header
    )
    
    
if __name__ == '__main__':

    main()
