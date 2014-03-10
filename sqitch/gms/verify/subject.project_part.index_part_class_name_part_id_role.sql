-- Verify subject.project_part.index_part_class_name_part_id_role

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'project_part_part_role_index';

ROLLBACK;
