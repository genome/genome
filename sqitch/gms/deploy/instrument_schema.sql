-- Deploy instrument_schema
-- requires: empty_db

BEGIN;

CREATE SCHEMA instrument AUTHORIZATION genome;

COMMIT;
