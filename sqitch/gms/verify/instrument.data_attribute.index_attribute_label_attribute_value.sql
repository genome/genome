-- Verify instrument.data_attribute.index_attribute_label_attribute_value

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_i_da_al_av';

ROLLBACK;
