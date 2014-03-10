-- Verify subject.project_part.index_label

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_s_pp_l';

ROLLBACK;
