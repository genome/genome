-- Verify timeline.analysis_project.pkey_id.sql

BEGIN;
SELECT 1/COUNT(pg_attribute.attname)
FROM pg_index, pg_class, pg_attribute
WHERE
  pg_class.oid = 'timeline.analysis_project'::regclass AND
  indrelid = pg_class.oid AND
  pg_attribute.attrelid = pg_class.oid AND
  pg_attribute.attnum = any(pg_index.indkey)
  AND indisprimary;
ROLLBACK;
