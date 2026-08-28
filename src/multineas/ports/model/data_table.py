from dataclasses import dataclass
import pandas as pd

@dataclass
class DataTable:
    data: pd.DataFrame
    source: str
    metadata: dict | None = None