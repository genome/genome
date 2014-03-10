-- Revert subject.subject_attribute.index_attribute_label

BEGIN;

DROP INDEX subject.subject_attribute_label_index;

COMMIT;
