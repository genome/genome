-- Deploy config_schema
-- requires: empty_db

BEGIN;

CREATE SCHEMA disk AUTHORIZATION genome;

COMMIT;
