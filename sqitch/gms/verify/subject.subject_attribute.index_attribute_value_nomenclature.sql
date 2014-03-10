-- Verify subject.subject_attribute.index_attribute_value_nomenclature

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_s_sa_av_n';

ROLLBACK;
