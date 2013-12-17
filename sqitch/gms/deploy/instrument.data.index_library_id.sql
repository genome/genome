-- Deploy instrument.data.library_id
-- requires: instrument_data

BEGIN;

CREATE INDEX instrument_data_library_id_index on instrument.data using btree (library_id);

COMMIT;
