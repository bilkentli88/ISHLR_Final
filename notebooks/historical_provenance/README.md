# Historical acquisition and analysis notebooks

These notebooks are retained to document the original iNaturalist data acquisition and validation process. They preserve the historical procedural source, but may require live network access or interactive execution and are **not** needed for the offline reproduction of archived results.

The canonical no-network workflow is `python scripts/reproduce_frozen_inaturalist.py`. Unlike an API rerun, it uses saved screening/event/trajectory tables that correspond to the published selected cohort. The underlying original API-response snapshot was not archived in full, and a new live API query may produce a different candidate set.
