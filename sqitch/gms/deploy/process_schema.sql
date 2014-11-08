-- Deploy process_schema
-- requires: empty_db

BEGIN;

CREATE SCHEMA process AUTHORIZATION genome;

COMMIT;
