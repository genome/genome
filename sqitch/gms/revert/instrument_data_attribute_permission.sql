-- Revert instrument_data_attribute_permission

BEGIN;

REVOKE ALL ON TABLE instrument.data_attribute FROM "gms-user";

COMMIT;
