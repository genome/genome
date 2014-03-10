-- Revert config_set_permission

BEGIN;

REVOKE ALL ON TABLE config.set FROM "gms-user";

COMMIT;
