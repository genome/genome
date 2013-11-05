-- Deploy subject.misc_attribute.entity_class_name_property_name_value
-- requires: subject_misc_attribute

BEGIN;

CREATE INDEX misc_attribute_entity_class_property_value_index on subject.misc_attribute using btree (entity_class_name, property_name, value);

COMMIT;
