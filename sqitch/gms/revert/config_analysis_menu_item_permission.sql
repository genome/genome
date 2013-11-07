-- Revert config_analysis_menu_item_permission

BEGIN;

REVOKE ALL ON TABLE config.analysis_menu_item FROM "gms-user";

COMMIT;
