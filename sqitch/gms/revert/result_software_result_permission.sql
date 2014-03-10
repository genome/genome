-- Revert result_software_result_permission

BEGIN;

REVOKE ALL ON TABLE result.software_result FROM "gms-user";

COMMIT;
