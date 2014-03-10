-- Deploy model_schema
-- requires: empty_db

BEGIN;

CREATE SCHEMA model AUTHORIZATION genome;

COMMIT;
