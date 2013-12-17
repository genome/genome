-- Verify instrument_data

BEGIN;

SELECT id, subclass_name, sequencing_platform, library_id,
    source_name, subset_name, run_name
FROM instrument.data
WHERE FALSE;

ROLLBACK;
