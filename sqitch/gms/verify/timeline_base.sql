-- Verify timeline_base

BEGIN;

SELECT id, created_by, updated_at, created_at, name, object_id, object_class_name, reason
FROM timeline.base
WHERE FALSE;

ROLLBACK;
