-- Verify instrument.data_attribute.index_attribute_label

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'instrument_data_attribute_label_index';

ROLLBACK;
