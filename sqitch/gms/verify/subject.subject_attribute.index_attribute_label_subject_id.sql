-- Verify subject.subject_attribute.index_attribute_label_subject_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_s_sa_al_si';

ROLLBACK;
