-- Deploy web_schema
-- requires: empty_db

BEGIN;

CREATE SCHEMA web AUTHORIZATION genome;

COMMIT;
