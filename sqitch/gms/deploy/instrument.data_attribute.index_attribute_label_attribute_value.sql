-- Deploy instrument.data_attribute.attribute_label_attribute_value
-- requires: instrument_data_attribute

BEGIN;

CREATE INDEX idx_i_da_al_av on instrument.data_attribute using btree (attribute_label, attribute_value);

COMMIT;
