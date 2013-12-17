-- Verify subject_misc_attribute

BEGIN;

SELECT entity_class_name, entity_id, property_name, value
FROM subject.misc_attribute
WHERE FALSE;

ROLLBACK;
