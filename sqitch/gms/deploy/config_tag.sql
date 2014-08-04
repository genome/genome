-- Deploy config_tag
-- requires: config_schema

BEGIN;
CREATE TABLE IF NOT EXISTS config.tag (
    id character varying(64) PRIMARY KEY,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    created_by character varying(255) NOT NULL,
    name text NOT NULL UNIQUE,
    description text
);

REVOKE ALL ON TABLE config.tag FROM PUBLIC;
GRANT ALL ON TABLE config.tag TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE config.tag TO "gms-user";
COMMIT;
