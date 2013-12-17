-- Deploy disk_schema_permissions
-- requires: disk_schema

BEGIN;

REVOKE ALL ON SCHEMA disk FROM PUBLIC;
REVOKE ALL ON SCHEMA disk FROM genome;
GRANT ALL ON SCHEMA disk TO genome;
GRANT USAGE ON SCHEMA disk TO "gms-user";

COMMIT;
