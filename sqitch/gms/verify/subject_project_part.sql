-- Verify subject_project_part

BEGIN;

SELECT id, project_id, part_class_name, part_id, label, role
FROM subject.project_part
WHERE FALSE;

ROLLBACK;
