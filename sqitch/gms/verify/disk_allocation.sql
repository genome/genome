-- Verify disk_allocation

BEGIN;

-- verify drop column preserved

SELECT 1/is_deleted
FROM (
    SELECT
        CASE WHEN count(*)=0 THEN 1
             ELSE 0
        END AS is_deleted
    FROM pg_catalog.pg_class c
    LEFT JOIN pg_catalog.pg_namespace n ON n.oid = c.relnamespace
    LEFT JOIN pg_catalog.pg_attribute a ON c.oid = a.attrelid
    WHERE c.relname ~ '^(allocation)$'
      AND n.nspname ~ '^(disk)$'
      AND a.attname ~ '^(preserved)$'
) q;

-- verify NOT NULL added to archive_after_time

SELECT 1/a.attnotnull::int
FROM pg_catalog.pg_class c
LEFT JOIN pg_catalog.pg_namespace n ON n.oid = c.relnamespace
LEFT JOIN pg_catalog.pg_attribute a ON c.oid = a.attrelid
WHERE c.relname ~ '^(allocation)$'
  AND n.nspname ~ '^(disk)$'
  AND a.attname ~ '^(archive_after_time)$'
;

-- verify NOT NULL added to status

SELECT 1/a.attnotnull::int
FROM pg_catalog.pg_class c
LEFT JOIN pg_catalog.pg_namespace n ON n.oid = c.relnamespace
LEFT JOIN pg_catalog.pg_attribute a ON c.oid = a.attrelid
WHERE c.relname ~ '^(allocation)$'
  AND n.nspname ~ '^(disk)$'
  AND a.attname ~ '^(status)$'
;




ROLLBACK;
