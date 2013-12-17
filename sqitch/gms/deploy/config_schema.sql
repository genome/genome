-- Deploy config_schema
-- requires: empty_db

BEGIN;

CREATE SCHEMA config AUTHORIZATION genome;

COMMIT;
