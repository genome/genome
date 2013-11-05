-- Verify subject.project_part.index_part_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_subject_project_part_part_id';

ROLLBACK;
