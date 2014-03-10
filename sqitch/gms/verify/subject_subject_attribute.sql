-- Verify subject_subject_attribute

BEGIN;

SELECT subject_id, attribute_label, attribute_value, nomenclature
FROM subject.subject_attribute
WHERE FALSE;

ROLLBACK;
