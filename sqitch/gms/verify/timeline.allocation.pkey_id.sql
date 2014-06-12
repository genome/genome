-- Verify timeline.allocation.pkey_id

BEGIN;
SELECT 1/COUNT(pg_attribute.attname)
FROM pg_index, pg_class, pg_attribute
WHERE
  pg_class.oid = 'timeline.allocation'::regclass AND
  indrelid = pg_class.oid AND
  pg_attribute.attrelid = pg_class.oid AND
  pg_attribute.attnum = any(pg_index.indkey)
  AND indisprimary;
ROLLBACK;
