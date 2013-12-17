-- Deploy instrument.data_attribute.attribute_label
-- requires: instrument_data_attribute

BEGIN;

CREATE INDEX instrument_data_attribute_label_index on instrument.data_attribute using btree (attribute_label);

COMMIT;
