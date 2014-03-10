-- Verify instrument_fragment_library

BEGIN;

SELECT library_id, full_name, sample_id, library_insert_size,
    original_insert_size, protocol, transcript_strand
FROM instrument.fragment_library
WHERE FALSE;

ROLLBACK;
