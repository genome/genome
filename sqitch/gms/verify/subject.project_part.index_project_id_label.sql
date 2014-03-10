-- Verify subject.project_part.index_project_id_label

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'project_part_project_label_index';

ROLLBACK;
