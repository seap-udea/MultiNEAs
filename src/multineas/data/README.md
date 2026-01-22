# Data directory for MultiNEAs
This directory contains data files used by the MultiNEAs package.

## Available Data

*   `nea_extended.json.gz`: Extended Near-Earth Asteroid dataset.

## Updating Data

To update the data files, run the update script from the project root:

```bash
./bin/data_update.sh
```

This script will download the latest version of `nea_extended.json.gz` from the [Minor Planet Center](https://minorplanetcenter.net/Extended_Files/nea_extended.json.gz).
