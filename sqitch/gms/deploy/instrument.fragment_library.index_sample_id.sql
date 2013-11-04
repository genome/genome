-- Deploy instrument.fragment_library.sample_id
-- requires: instrument_fragment_library

BEGIN;

CREATE INDEX fragment_library_sample_id_index on instrument.fragment_library using btree (sample_id);

COMMIT;
