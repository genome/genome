-- Verify subject.misc_attribute.index_entity_class_name_property_name_value

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'misc_attribute_entity_class_property_value_index';

ROLLBACK;
