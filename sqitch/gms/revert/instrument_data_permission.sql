-- Revert instrument_data_permission

BEGIN;

REVOKE ALL ON TABLE instrument.data FROM "gms-user";

COMMIT;
