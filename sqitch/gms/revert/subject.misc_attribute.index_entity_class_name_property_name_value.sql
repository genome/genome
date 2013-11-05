-- Revert subject.misc_attribute.index_entity_class_name_property_name_value

BEGIN;

DROP INDEX subject.misc_attribute_entity_class_property_value_index;

COMMIT;
