-- Revert result_input_permission

BEGIN;

REVOKE ALL ON TABLE result.input FROM "gms-user";

COMMIT;
