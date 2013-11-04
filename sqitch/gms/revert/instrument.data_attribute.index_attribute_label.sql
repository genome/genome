-- Revert instrument.data_attribute.index_attribute_label

BEGIN;

DROP INDEX instrument.instrument_data_attribute_label_index;

COMMIT;
