-- Verify subject.project.index_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'project_name_index';

ROLLBACK;
