-- Verify subject.project_part.index_project_id_role

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'project_part_project_role_index';

ROLLBACK;
