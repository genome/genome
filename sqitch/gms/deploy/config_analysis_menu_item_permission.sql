-- Deploy config_analysis_menu_item_permission
-- requires: config_analysis_menu_item

BEGIN;

REVOKE ALL ON TABLE config.analysis_menu_item FROM PUBLIC;
REVOKE ALL ON TABLE config.analysis_menu_item FROM genome;
GRANT ALL ON TABLE config.analysis_menu_item TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE config.analysis_menu_item TO "gms-user";

COMMIT;
