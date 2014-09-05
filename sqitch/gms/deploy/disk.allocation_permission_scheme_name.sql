-- Deploy disk.allocation_permission_scheme_name
-- requires: disk.allocation_permission_scheme

BEGIN;

ALTER TABLE disk.allocation_permission_scheme ADD COLUMN name character varying(128);

COMMIT;
