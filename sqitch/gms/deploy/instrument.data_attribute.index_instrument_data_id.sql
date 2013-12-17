-- Deploy instrument.data_attribute.instrument_data_id
-- requires: instrument_data_attribute

BEGIN;

CREATE INDEX instrument_data_id_index on instrument.data_attribute using btree (instrument_data_id);

COMMIT;
