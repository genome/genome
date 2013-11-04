-- Deploy model.event.instrument_data_id
-- requires: model_event

BEGIN;

CREATE INDEX event_inst_data_index on model.event using btree (instrument_data_id);

COMMIT;
