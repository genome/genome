-- Revert model.event.index_instrument_data_id

BEGIN;

DROP INDEX model.event_inst_data_index;

COMMIT;
