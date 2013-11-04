-- Revert instrument.data_attribute.index_attribute_label_attribute_value

BEGIN;

DROP INDEX instrument.idx_i_da_al_av;

COMMIT;
