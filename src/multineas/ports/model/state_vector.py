from dataclasses import dataclass
from typing import List

@dataclass
class StateVector:
    body: str
    epoch: str
    position: List[float]   # [x, y, z] en km
    velocity: List[float]   # [vx, vy, vz] en km/s