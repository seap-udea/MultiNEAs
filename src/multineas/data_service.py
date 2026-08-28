from multineas.ports.model.data_table import DataTable
from multineas.ports.adapter.cneos_adapter import CNEOSAdapter
from multineas.ports.adapter.sbdb_adapter import SBDBAdapter

def get_data(source: str, params: dict) -> DataTable:
    INPUT_ADAPTERS = {
        'CNEOS':CNEOSAdapter, 
        'SBDB':SBDBAdapter
    }

    try:
        adapter = INPUT_ADAPTERS[source]()
    except KeyError:
        raise ValueError(f"Unsupported source: {source}, available sources: {INPUT_ADAPTERS.keys()}")

    return adapter.get_data(params)