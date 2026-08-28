from abc import ABC, abstractmethod
from multineas.ports.model.data_table import DataTable

class DatabaseConnectorPort(ABC):
    @abstractmethod
    def get_data(self, params: dict) -> DataTable:
        pass