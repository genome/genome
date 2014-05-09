-- Revert remove_unused_apconfig_tables

BEGIN;

CREATE TABLE IF NOT EXISTS config.set (
  id character varying(64) NOT NULL,
  created_at timestamp(6) without time zone NOT NULL,
  updated_at timestamp(6) without time zone NOT NULL,
  allocation_id character varying(64) NOT NULL,
  CONSTRAINT genome_config_set_pk PRIMARY KEY (id),
  CONSTRAINT config_set_allocation_fk FOREIGN KEY (allocation_id) REFERENCES disk.allocation(id)
);
REVOKE ALL ON TABLE config.set FROM PUBLIC;
REVOKE ALL ON TABLE config.set FROM genome;
GRANT ALL ON TABLE config.set TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE config.set TO "gms-user";

CREATE TABLE IF NOT EXISTS config.analysis_menu_item (
  id character varying(64) NOT NULL,
  created_at timestamp(6) without time zone NOT NULL,
  updated_at timestamp(6) without time zone NOT NULL,
  name character varying(255) NOT NULL,
  configuration_set_id character varying(64) NOT NULL,
  CONSTRAINT analysis_menu_item_pk PRIMARY KEY (id),
  CONSTRAINT analysis_menu_config_set_fk FOREIGN KEY (configuration_set_id) REFERENCES config.set(id)
);
REVOKE ALL ON TABLE config.analysis_menu_item FROM PUBLIC;
REVOKE ALL ON TABLE config.analysis_menu_item FROM genome;
GRANT ALL ON TABLE config.analysis_menu_item TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE config.set TO "gms-user";

ALTER TABLE config.analysis_project ADD COLUMN analysis_menu_item_id character varying(64) REFERENCES config.analysis_menu_item(id);
ALTER TABLE config.analysis_project ADD COLUMN configuration_set_id character varying(64) REFERENCES config.set(id);

COMMIT;
