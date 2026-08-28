from multineas.ports.model.data_table import DataTable
from multineas.ports.database_conector_port import DatabaseConnectorPort
import requests
import pandas as pd

class CNEOSAdapter(DatabaseConnectorPort): 

    BASE_URL = "https://ssd-api.jpl.nasa.gov/fireball.api"

    valid_params = [
        "date-min",
        "date-max",
        "energy-min",
        "energy-max",
        "impact-e-min",
        "impact-e-max",
        "alt-min",
        "alt-max",
        "req-loc",
        "req-alt",
        "req-vel-comp",
        "vel-comp",
        "sort",
        "limit",
    ]

    def get_data(self, params: dict) -> DataTable:
        invalid = set(params.keys()) - set(self.valid_params)
        if invalid:
            raise ValueError(
                f"Invalid parameter(s): {sorted(invalid)}. "
                f"Valid parameters for CNEOS: {sorted(self.valid_params)}"
            )

        r = requests.get(self.BASE_URL, params=params)
        response = r.json()

        return DataTable(
            data=pd.DataFrame(response["data"], columns=response["fields"]), 
            source="CNEOS",
            metadata=params
        )