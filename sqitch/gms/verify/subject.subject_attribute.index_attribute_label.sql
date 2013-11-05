-- Verify subject.subject_attribute.index_attribute_label

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'subject_attribute_label_index';

ROLLBACK;
